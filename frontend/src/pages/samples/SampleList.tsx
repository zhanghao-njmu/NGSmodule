/**
 * Sample List Page - Sample management with CSV import
 */
import React, { useEffect, useState } from 'react'
import {
  Card,
  Button,
  Table,
  Space,
  Select,
  Upload,
  message,
  Tag,
  Tooltip,
} from 'antd'
import {
  PlusOutlined,
  UploadOutlined,
  ExperimentOutlined,
  FileTextOutlined,
} from '@ant-design/icons'
import type { ColumnsType } from 'antd/es/table'
import { useSampleStore } from '../../store/sampleStore'
import { useProjectStore } from '../../store/projectStore'
import type { Sample } from '../../types/sample'
import dayjs from 'dayjs'

const { Option } = Select

export const SampleList: React.FC = () => {
  const { samples, loading, fetchSamples, importFromCSV, deleteSample } = useSampleStore()
  const { projects, fetchProjects } = useProjectStore()
  const [selectedProject, setSelectedProject] = useState<string>('')

  useEffect(() => {
    fetchProjects()
  }, [])

  useEffect(() => {
    if (selectedProject) {
      fetchSamples({ project_id: selectedProject })
    }
  }, [selectedProject])

  const handleCSVImport = (file: File) => {
    if (selectedProject) {
      importFromCSV(selectedProject, file)
    } else {
      message.warning('Please select a project first')
    }
    return false // Prevent auto upload
  }

  const columns: ColumnsType<Sample> = [
    {
      title: 'Sample Name',
      dataIndex: 'sample_name',
      key: 'sample_name',
      render: (name) => (
        <Space>
          <ExperimentOutlined style={{ color: '#52c41a' }} />
          <span style={{ fontWeight: 500 }}>{name}</span>
        </Space>
      ),
    },
    {
      title: 'Type',
      dataIndex: 'sample_type',
      key: 'sample_type',
      width: 150,
      render: (type) => type && <Tag color="blue">{type}</Tag>,
    },
    {
      title: 'Files',
      dataIndex: 'file_count',
      key: 'file_count',
      width: 100,
      render: (count) => count || 0,
    },
    {
      title: 'Created',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 150,
      render: (date) => dayjs(date).format('YYYY-MM-DD HH:mm'),
    },
  ]

  return (
    <div>
      <Card style={{ marginBottom: 16 }}>
        <Space style={{ width: '100%', justifyContent: 'space-between' }}>
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
          <Space>
            <Upload beforeUpload={handleCSVImport} accept=".csv" showUploadList={false}>
              <Button icon={<UploadOutlined />}>Import CSV</Button>
            </Upload>
            <Button type="primary" icon={<PlusOutlined />} disabled={!selectedProject}>
              New Sample
            </Button>
          </Space>
        </Space>
      </Card>

      <Card>
        <Table
          columns={columns}
          dataSource={samples}
          rowKey="id"
          loading={loading}
          pagination={{ pageSize: 20 }}
        />
      </Card>
    </div>
  )
}

export default SampleList
