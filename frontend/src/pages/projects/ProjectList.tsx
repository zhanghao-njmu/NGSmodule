/**
 * Project List Page
 */
import React from 'react'
import { Card, Button, Empty, Typography } from 'antd'
import { PlusOutlined, FolderOutlined } from '@ant-design/icons'

const { Title } = Typography

export const ProjectList: React.FC = () => {
  return (
    <div>
      <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: 24 }}>
        <Title level={2} style={{ margin: 0 }}>
          Projects
        </Title>
        <Button type="primary" icon={<PlusOutlined />} size="large">
          New Project
        </Button>
      </div>

      <Card>
        <Empty
          image={<FolderOutlined style={{ fontSize: 64, color: '#ccc' }} />}
          description="No projects yet. Create your first project to get started."
        >
          <Button type="primary" icon={<PlusOutlined />}>
            Create Project
          </Button>
        </Empty>
      </Card>
    </div>
  )
}

export default ProjectList
