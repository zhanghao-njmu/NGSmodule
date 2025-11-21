/**
 * File List Page - File upload and management
 */
import React, { useEffect, useState } from 'react'
import {
  Card,
  Button,
  Progress,
  Select,
  Tag,
} from 'antd'
import {
  UploadOutlined,
  DownloadOutlined,
  FileOutlined,
} from '@ant-design/icons'
import type { ColumnsType } from 'antd/es/table'
import { useFileStore } from '../../store/fileStore'
import { useProjectStore } from '../../store/projectStore'
import { PageHeader, DataTable } from '../../components/common'
import type { FileItem } from '../../types/file'
import dayjs from 'dayjs'

const { Option } = Select

export const FileList: React.FC = () => {
  const { files, loading, uploadProgress, fetchFiles, downloadFile } = useFileStore()
  const { projects, fetchProjects } = useProjectStore()
  const [selectedProject, setSelectedProject] = useState<string>('')

  useEffect(() => {
    fetchProjects()
  }, [])

  useEffect(() => {
    if (selectedProject) {
      fetchFiles({ project_id: selectedProject })
    }
  }, [selectedProject])

  const formatFileSize = (bytes: number) => {
    if (bytes === 0) return '0 B'
    const k = 1024
    const sizes = ['B', 'KB', 'MB', 'GB', 'TB']
    const i = Math.floor(Math.log(bytes) / Math.log(k))
    return `${(bytes / Math.pow(k, i)).toFixed(2)} ${sizes[i]}`
  }

  const columns: ColumnsType<FileItem> = [
    {
      title: 'Filename',
      dataIndex: 'filename',
      key: 'filename',
      render: (name) => (
        <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
          <FileOutlined style={{ color: 'var(--color-primary)' }} />
          <span>{name}</span>
        </div>
      ),
    },
    {
      title: 'Type',
      dataIndex: 'file_type',
      key: 'file_type',
      width: 120,
      render: (type) => type && <Tag>{type}</Tag>,
    },
    {
      title: 'Size',
      dataIndex: 'file_size',
      key: 'file_size',
      width: 120,
      render: (size) => formatFileSize(size),
    },
    {
      title: 'Uploaded',
      dataIndex: 'uploaded_at',
      key: 'uploaded_at',
      width: 150,
      render: (date) => dayjs(date).format('YYYY-MM-DD HH:mm'),
    },
    {
      title: 'Actions',
      key: 'actions',
      width: 100,
      render: (_, record) => (
        <Button
          icon={<DownloadOutlined />}
          size="small"
          onClick={() => downloadFile(record.id, record.filename)}
        >
          Download
        </Button>
      ),
    },
  ]

  return (
    <div>
      <PageHeader
        left={
          <Select
            placeholder="Select Project"
            style={{ width: 300 }}
            value={selectedProject || undefined}
            onChange={setSelectedProject}
          >
            {projects.map((p) => (
              <Option key={p.id} value={p.id}>
                {p.name}
              </Option>
            ))}
          </Select>
        }
        right={
          <Button type="primary" icon={<UploadOutlined />} disabled={!selectedProject}>
            Upload Files
          </Button>
        }
      />

      {uploadProgress > 0 && uploadProgress < 100 && (
        <Card style={{ marginBottom: 16 }}>
          <Progress percent={uploadProgress} />
        </Card>
      )}

      <DataTable
        columns={columns}
        dataSource={files}
        rowKey="id"
        loading={loading}
        pagination={{ pageSize: 20 }}
        emptyText="No Files"
        emptyDescription="Select a project and upload files to get started"
      />
    </div>
  )
}

export default FileList
