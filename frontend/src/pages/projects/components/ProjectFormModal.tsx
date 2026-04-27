/**
 * Project Form Modal — TanStack Query mutation version.
 */
import type React from 'react'
import { useEffect } from 'react'
import { Modal, Form, Input, Select } from 'antd'

import { useCreateProject, useUpdateProject } from '@/hooks/queries'
import type { Project } from '@/types/project'

const { TextArea } = Input
const { Option } = Select

interface ProjectFormModalProps {
  open: boolean
  project: Project | null
  onClose: () => void
  onSuccess: () => void
}

export const ProjectFormModal: React.FC<ProjectFormModalProps> = ({
  open,
  project,
  onClose,
  onSuccess,
}) => {
  const [form] = Form.useForm()
  const createMutation = useCreateProject()
  const updateMutation = useUpdateProject()
  const loading = createMutation.isPending || updateMutation.isPending

  useEffect(() => {
    if (open) {
      if (project) {
        form.setFieldsValue({
          name: project.name,
          description: project.description,
          status: project.status,
        })
      } else {
        form.resetFields()
      }
    }
  }, [open, project, form])

  const handleSubmit = async () => {
    try {
      const values = await form.validateFields()

      if (project) {
        await updateMutation.mutateAsync({ id: project.id, data: values })
      } else {
        await createMutation.mutateAsync(values)
      }

      onSuccess()
    } catch (error) {
      console.error('Form submission error:', error)
    }
  }

  return (
    <Modal
      title={project ? 'Edit Project' : 'Create New Project'}
      open={open}
      onCancel={onClose}
      onOk={handleSubmit}
      confirmLoading={loading}
      width={600}
      destroyOnClose
    >
      <Form form={form} layout="vertical" initialValues={{ status: 'active' }}>
        <Form.Item
          name="name"
          label="Project Name"
          rules={[
            { required: true, message: 'Please enter project name' },
            { min: 2, message: 'Name must be at least 2 characters' },
            { max: 100, message: 'Name must not exceed 100 characters' },
          ]}
        >
          <Input placeholder="Enter project name" />
        </Form.Item>

        <Form.Item
          name="description"
          label="Description"
          rules={[{ max: 500, message: 'Description must not exceed 500 characters' }]}
        >
          <TextArea
            rows={4}
            placeholder="Enter project description (optional)"
            showCount
            maxLength={500}
          />
        </Form.Item>

        {project && (
          <Form.Item
            name="status"
            label="Status"
            rules={[{ required: true, message: 'Please select status' }]}
          >
            <Select>
              <Option value="active">Active</Option>
              <Option value="completed">Completed</Option>
              <Option value="archived">Archived</Option>
            </Select>
          </Form.Item>
        )}
      </Form>
    </Modal>
  )
}

export default ProjectFormModal
